namespace PolyploidQtlSeqCore.VariantCall
{
    /// <summary>
    /// SnpEff設定
    /// </summary>
    internal class SnpEffSettings
    {
        /// <summary>
        /// SnpEff設定を作成する。
        /// </summary>
        /// <param name="settingValue">オプションの値</param>
        public SnpEffSettings(ISnpEffSettingValue settingValue)
        {
            MaxHeap = new SnpEffMaxHeap(settingValue.SnpEffMaxHeap);
            ConfigFile = new SnpEffConfigFile(settingValue.SnpEffConfigFile);
            Database = new SnpEffDatabase(settingValue.SnpEffDatabaseName);
        }


        /// <summary>
        /// MaxHeapサイズを取得する。
        /// </summary>
        public SnpEffMaxHeap MaxHeap { get; }

        /// <summary>
        /// SnpEff.configファイルを取得する。
        /// </summary>
        public SnpEffConfigFile ConfigFile { get; }

        /// <summary>
        /// SnpEffデータベースを取得する。
        /// </summary>
        public SnpEffDatabase Database { get; }

        /// <summary>
        /// SnpEffが実行可能かどうかを取得する。
        /// </summary>
        public bool CanSneEff => !string.IsNullOrEmpty(Database.Value);
    }
}
