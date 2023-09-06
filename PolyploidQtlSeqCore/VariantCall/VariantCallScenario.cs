using Kurukuru;
using PolyploidQtlSeqCore.IO;
using PolyploidQtlSeqCore.Mapping;
using PolyploidQtlSeqCore.QtlAnalysis;
using PolyploidQtlSeqCore.QtlAnalysis.Chr;

namespace PolyploidQtlSeqCore.VariantCall
{
    /// <summary>
    /// 変異検出シナリオ
    /// </summary>
    internal class VariantCallScenario
    {
        private const string SNPEFF_VCF_FILENAME = "polyQtlseq.snpEff.vcf.gz";

        private readonly AnalysisChrSettings _analysisChrSettings;
        private readonly BcfToolsVariantCallSettings _variantCallSettings;
        private readonly SnpEffSettings _snpEffSettings;


        /// <summary>
        /// 変異検出シナリオインスタンスを作成する。
        /// </summary>
        /// <param name="settings">変異検出シナリオ設定</param>
        public VariantCallScenario(VariantCallScenarioSettings settings)
        {
            _analysisChrSettings = settings.AnalysisChrSettings;
            _variantCallSettings= settings.BcfToolsVariantCallSettings;
            _snpEffSettings = settings.SnpEffSettings;
        }

        /// <summary>
        /// 指定染色体の変異検出+SnpEffを行う。
        /// </summary>
        /// <param name="bamFiles">全サンプルのBAMファイル</param>
        /// <returns>VCFファイル</returns>
        public async ValueTask<VcfFile> CallAsync(AllSampleBamFiles bamFiles)
        {
            Log.Clear();
            var logFilePath = _variantCallSettings.OutputDirectory.CreateFilePath("VariantCall.Log.txt");
            var analysisChrs = bamFiles.GetAnalysisChrs(_analysisChrSettings);

            try
            {
                _variantCallSettings.OutputDirectory.Create();
                var qtlseqVcf = await ParallelCallAsync(bamFiles, analysisChrs);
                if (!_snpEffSettings.CanSneEff) return qtlseqVcf;

                var qtlseqSnpEffVcf = await SnpEffAsync(qtlseqVcf, _variantCallSettings.OutputDirectory);
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


        private async ValueTask<VcfFile> ParallelCallAsync(AllSampleBamFiles bamFiles, Chromosome[] chromosomes)
        {
            return await Spinner.StartAsync("Variant call ...", async spinner =>
            {
                try
                {
                    var variantCaller = new ParallelVariantCall(_variantCallSettings);
                    var qtlseqVcf = await variantCaller.CallAsync(bamFiles, chromosomes);
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
                    var snpEffPipeline = new SnpEffPipeline(_snpEffSettings);
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
