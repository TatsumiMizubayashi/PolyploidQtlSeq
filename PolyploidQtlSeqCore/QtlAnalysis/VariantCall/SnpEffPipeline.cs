namespace PolyploidQtlSeqCore.QtlAnalysis.VariantCall
{
    /// <summary>
    /// SnpEff, VCF圧縮, indexファイル作成パイプライン
    /// </summary>
    internal class SnpEffPipeline
    {
        private const string GzExtension = ".gz";
        private readonly SnpEffOption _option;

        /// <summary>
        /// SnpEffパイプラインを作成する。
        /// </summary>
        /// <param name="option">オプション</param>
        public SnpEffPipeline(SnpEffOption option)
        {
            _option = option;
        }

        /// <summary>
        /// SnpEffパイプラインを実行する。
        /// </summary>
        /// <param name="inputVcf">入力VCF</param>
        /// <param name="outputVcfPath">出力VCFファイルPath</param>
        /// <returns>SnpEff済みVCFファイル</returns>
        public async ValueTask<VcfFile> RunAsync(VcfFile inputVcf, string outputVcfPath)
        {
            var snpEff = new SnpEff(_option);
            var snpEffOutputVcfPath = GetSnpEffVcfPath(outputVcfPath);
            var snpEffVcf = await snpEff.RunAsync(inputVcf, snpEffOutputVcfPath);
            var snpEffCompressedVcf = await snpEffVcf.CompressionAsync();
            await snpEffCompressedVcf.CreateIndexFileAsync();

            return snpEffCompressedVcf;
        }

        private static string GetSnpEffVcfPath(string outputVcfPath)
        {
            var extension = Path.GetExtension(outputVcfPath);
            if (extension != GzExtension) return outputVcfPath;

            var saveDirPath = Path.GetDirectoryName(outputVcfPath);
            if (string.IsNullOrEmpty(saveDirPath)) throw new ArgumentException(null, nameof(outputVcfPath));

            return Path.Combine(saveDirPath, Path.GetFileNameWithoutExtension(outputVcfPath));
        }
    }
}
