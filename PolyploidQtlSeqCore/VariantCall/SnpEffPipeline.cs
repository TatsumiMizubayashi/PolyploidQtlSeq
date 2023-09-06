namespace PolyploidQtlSeqCore.VariantCall
{
    /// <summary>
    /// SnpEff, VCF圧縮, indexファイル作成パイプライン
    /// </summary>
    internal class SnpEffPipeline
    {
        private const string GzExtension = ".gz";
        private readonly SnpEffSettings _settings;

        /// <summary>
        /// SnpEffパイプラインを作成する。
        /// </summary>
        /// <param name="settings">設定</param>
        public SnpEffPipeline(SnpEffSettings settings)
        {
            _settings = settings;
        }

        /// <summary>
        /// SnpEffパイプラインを実行する。
        /// </summary>
        /// <param name="inputVcf">入力VCF</param>
        /// <param name="outputVcfPath">出力VCFファイルPath</param>
        /// <returns>SnpEff済みVCFファイル</returns>
        public async ValueTask<VcfFile> RunAsync(VcfFile inputVcf, string outputVcfPath)
        {
            var snpEff = new SnpEff(_settings);
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
