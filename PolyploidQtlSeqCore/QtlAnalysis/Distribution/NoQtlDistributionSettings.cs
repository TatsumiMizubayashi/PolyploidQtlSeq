namespace PolyploidQtlSeqCore.QtlAnalysis.Distribution
{
    /// <summary>
    /// QTLなし分布作成設定
    /// </summary>
    internal class NoQtlDistributionSettings
    {
        /// <summary>
        /// QTLなし分布作成設定を作成する。
        /// </summary>
        /// <param name="settingValue">設定値値</param>
        /// <param name="parameterDictionary">パラメータファイルの中身</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public NoQtlDistributionSettings(INoQtlDistributionSettingValue settingValue)
        {
            Ploidy = new Ploidy(settingValue.Ploidy);
            Parent2PlexNumber = new Parent2PlexNumber(settingValue.Parent2PlexNumber);
            Bulk1Number = new Bulk1Number(settingValue.Bulk1Number);
            Bulk2Number = new Bulk2Number(settingValue.Bulk2Number);
            ReplicatesNumber = new ReplicatesNumber(settingValue.ReplicatesNumber);
        }

        /// <summary>
        /// 倍数性を取得する。
        /// </summary>
        public Ploidy Ploidy { get; }

        /// <summary>
        /// 親2のPlex数を取得する。
        /// </summary>
        public Parent2PlexNumber Parent2PlexNumber { get; }

        /// <summary>
        /// Bulk1個体数を取得する。
        /// </summary>
        public Bulk1Number Bulk1Number { get; }

        /// <summary>
        /// Bulk2個体数を取得する。
        /// </summary>
        public Bulk2Number Bulk2Number { get; }

        /// <summary>
        /// 分布作成時の試行回数を取得する。
        /// </summary>
        public ReplicatesNumber ReplicatesNumber { get; }
    }
}
