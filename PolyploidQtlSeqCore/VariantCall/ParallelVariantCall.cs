using PolyploidQtlSeqCore.Mapping;
using PolyploidQtlSeqCore.QtlAnalysis.Chr;

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

        /// <summary>
        /// 並列変異検出を作成する。
        /// </summary>
        /// <param name="settings">設定</param>
        /// <param name="thread">スレッド数</param>
        public ParallelVariantCall(BcfToolsVariantCallSettings settings)
        {
            _settings = settings;
        }

        /// <summary>
        /// 変異検出を行う。
        /// </summary>
        /// <param name="bamFiles">BAMファイル</param>
        /// <param name="chromosomes">染色体</param>
        /// <returns>VCFファイル</returns>
        public async ValueTask<VcfFile> CallAsync(AllSampleBamFiles bamFiles, Chromosome[] chromosomes)
        {
            var pOption = new ParallelOptions()
            {
                MaxDegreeOfParallelism = _settings.ThreadNumber.Value
            };

            var chrVcfList = new List<OneChromosomeVcfFile>();
            await Parallel.ForEachAsync(chromosomes, pOption, async (chr, cancelToken) =>
            {
                var variantCallPipeline = new BcftoolsVariantCallPipeline(_settings);
                var chrVcfFile = await variantCallPipeline.CallAsync(bamFiles, chr);
                lock (_syncObj) chrVcfList.Add(chrVcfFile);
            });

            var mergeVcfFilePath = _settings.OutputDirectory.CreateFilePath(VCF_FILENAME);
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
