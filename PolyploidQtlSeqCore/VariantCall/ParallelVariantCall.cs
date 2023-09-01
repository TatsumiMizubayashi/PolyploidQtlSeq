using PolyploidQtlSeqCore.Mapping;
using PolyploidQtlSeqCore.QtlAnalysis;
using PolyploidQtlSeqCore.QtlAnalysis.Chr;
using PolyploidQtlSeqCore.Share;

namespace PolyploidQtlSeqCore.VariantCall
{
    /// <summary>
    /// 並列変異検出
    /// </summary>
    internal class ParallelVariantCall
    {
        private const string VCF_FILENAME = "polyQtlseq.vcf.gz";

        private readonly object _syncObj = new();
        private readonly BcfToolsVariantCallSettings _settings;
        private readonly ThreadNumber _threadNumber;

        /// <summary>
        /// 並列変異検出を作成する。
        /// </summary>
        /// <param name="option">オプション</param>
        /// <param name="thread">スレッド数</param>
        public ParallelVariantCall(BcfToolsVariantCallSettings option, ThreadNumber thread)
        {
            _settings = option;
            _threadNumber = thread;
        }

        /// <summary>
        /// 変異検出を行う。
        /// </summary>
        /// <param name="bamFiles">BAMファイル</param>
        /// <param name="outputDir">出力ディレクトリ</param>
        /// <param name="chromosomes">染色体</param>
        /// <returns>VCFファイル</returns>
        public async ValueTask<VcfFile> CallAsync(AllSampleBamFiles bamFiles, OutputDirectory outputDir, Chromosome[] chromosomes)
        {
            var pOption = new ParallelOptions()
            {
                MaxDegreeOfParallelism = _threadNumber.Value
            };

            var chrVcfList = new List<OneChromosomeVcfFile>();
            await Parallel.ForEachAsync(chromosomes, pOption, async (chr, cancelToken) =>
            {
                var variantCallPipeline = new BcftoolsVariantCallPipeline(_settings);
                var chrVcfFile = await variantCallPipeline.CallAsync(bamFiles, outputDir, chr);
                lock (_syncObj) chrVcfList.Add(chrVcfFile);
            });

            var mergeVcfFilePath = outputDir.CreateFilePath(VCF_FILENAME);
            var vcfFile = await BcftoolsConcat.RunAsync(mergeVcfFilePath, chrVcfList);
            await vcfFile.CreateIndexFileAsync();

            foreach (var chrVcf in chrVcfList)
            {
                chrVcf.Delete();
            }

            return vcfFile;
        }
    }
}
