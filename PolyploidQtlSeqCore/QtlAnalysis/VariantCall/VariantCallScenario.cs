using Kurukuru;
using PolyploidQtlSeqCore.IO;
using PolyploidQtlSeqCore.QtlAnalysis.Chr;
using PolyploidQtlSeqCore.QtlAnalysis.Mapping;

namespace PolyploidQtlSeqCore.QtlAnalysis.VariantCall
{
    /// <summary>
    /// 変異検出シナリオ
    /// </summary>
    internal class VariantCallScenario
    {
        private const string SNPEFF_VCF_FILENAME = "polyQtlseq.snpEff.vcf.gz";

        private readonly BcfToolsVariantCallOption _bcfToolsVariantCallOption;
        private readonly SnpEffOption _snpEffOption;
        private readonly ReferenceSequence _refSeq;
        private readonly ThreadNumber _thread;

        /// <summary>
        /// 変異検出シナリオを作成する。
        /// </summary>
        /// <param name="bcftoolsOption">bcftoolsオプション</param>
        /// <param name="snpEffOption">snpEffオプション</param>
        /// <param name="refSeq">リファレンスシークエンス</param>
        /// <param name="thread">使用するスレッド数</param>
        public VariantCallScenario(BcfToolsVariantCallOption bcftoolsOption, SnpEffOption snpEffOption,
            ReferenceSequence refSeq, ThreadNumber thread)
        {
            _bcfToolsVariantCallOption = bcftoolsOption;
            _snpEffOption = snpEffOption;
            _refSeq = refSeq;
            _thread = thread;
        }

        /// <summary>
        /// 指定染色体の変異検出+SnpEffを行う。
        /// </summary>
        /// <param name="bamFiles">全サンプルのBAMファイル</param>
        /// <param name="outputDir">出力ディレクトリ</param>
        /// <param name="chromosomes">解析する染色体</param>
        /// <returns>VCFファイル</returns>
        public async ValueTask<VcfFile> CallAsync(AllSampleBamFiles bamFiles, OutputDirectory outputDir, Chromosome[] chromosomes)
        {
            Log.Clear();
            var logFilePath = outputDir.CreateFilePath("VariantCall.Log.txt");

            try
            {
                var qtlseqVcf = await ParallelCallAsync(bamFiles, outputDir, chromosomes);
                if (!_snpEffOption.CanSneEff) return qtlseqVcf;

                var qtlseqSnpEffVcf = await SnpEffAsync(qtlseqVcf, outputDir);
                qtlseqVcf.Delete();

                return qtlseqSnpEffVcf;
            }
            catch
            {
                throw;
            }
            finally
            {
                Log.Save(logFilePath);
            }
        }

        private async ValueTask<VcfFile> ParallelCallAsync(AllSampleBamFiles bamFiles, OutputDirectory outputDir, Chromosome[] chromosomes)
        {
            return await Spinner.StartAsync("Variant call ...", async spinner =>
            {
                try
                {
                    var variantCaller = new ParallelVariantCall(_bcfToolsVariantCallOption, _refSeq, _thread);
                    var qtlseqVcf = await variantCaller.CallAsync(bamFiles, outputDir, chromosomes);
                    spinner.Succeed("Variant call completed");

                    return qtlseqVcf;
                }
                catch
                {
                    spinner.Fail("Variant call error");

                    throw;
                }
            });
        }

        private async ValueTask<VcfFile> SnpEffAsync(VcfFile inputVcf, OutputDirectory outputDir)
        {
            return await Spinner.StartAsync("SnpEff ...", async spinner =>
            {
                try
                {
                    var snpEffPipeline = new SnpEffPipeline(_snpEffOption);
                    var snpEffVcfPath = outputDir.CreateFilePath(SNPEFF_VCF_FILENAME);
                    var qtlseqSnpEffVcf = await snpEffPipeline.RunAsync(inputVcf, snpEffVcfPath);
                    spinner.Succeed("SnpEff completed");

                    return qtlseqSnpEffVcf;
                }
                catch
                {
                    spinner.Fail("SnpEff error");

                    throw;
                }
            });
        }
    }
}
