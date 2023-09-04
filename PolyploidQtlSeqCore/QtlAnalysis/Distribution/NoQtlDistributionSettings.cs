namespace PolyploidQtlSeqCore.QtlAnalysis.Distribution
{
    /// <summary>
    /// QTLなし分布作成設定
    /// </summary>
    internal class NoQtlDistributionSettings
    {
        private static readonly IReadOnlyDictionary<string, string> _toLongNameDictionary = new Dictionary<string, string>()
        {
            [Ploidy.SHORT_NAME] = Ploidy.LONG_NAME,
            [Ploidy.LONG_NAME] = Ploidy.LONG_NAME,

            [Parent2PlexNumber.SHORT_NAME] = Parent2PlexNumber.LONG_NAME,
            [Parent2PlexNumber.LONG_NAME] = Parent2PlexNumber.LONG_NAME,

            [Bulk1Number.SHORT_NAME] = Bulk1Number.LONG_NAME,
            [Bulk1Number.LONG_NAME] = Bulk1Number.LONG_NAME,

            [Bulk2Number.SHORT_NAME] = Bulk2Number.LONG_NAME,
            [Bulk2Number.LONG_NAME] = Bulk2Number.LONG_NAME,

            [ReplicatesNumber.SHORT_NAME] = ReplicatesNumber.LONG_NAME,
            [ReplicatesNumber.LONG_NAME] = ReplicatesNumber.LONG_NAME
        };

        /// <summary>
        /// LongName変換項目を辞書に追加する。
        /// </summary>
        /// <param name="dictionary">LongName変換辞書</param>
        public static void AddLongNameKeyValuePair(Dictionary<string, string> dictionary)
        {
            foreach (var keyValuePair in _toLongNameDictionary)
            {
                dictionary.Add(keyValuePair.Key, keyValuePair.Value);
            }
        }

        /// <summary>
        /// QTLなし分布作成設定を作成する。
        /// </summary>
        /// <param name="settingValue">設定値値</param>
        /// <param name="parameterDictionary">パラメータファイルの中身</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        [Obsolete("削除予定")]
        public NoQtlDistributionSettings(INoQtlDistributionSettingValue settingValue, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            Ploidy = new Ploidy(settingValue.Ploidy, parameterDictionary, userOptionDictionary);
            Parent2PlexNumber = new Parent2PlexNumber(settingValue.Parent2PlexNumber, parameterDictionary, userOptionDictionary);
            Bulk1Number = new Bulk1Number(settingValue.Bulk1Number, parameterDictionary, userOptionDictionary);
            Bulk2Number = new Bulk2Number(settingValue.Bulk2Number, parameterDictionary, userOptionDictionary);
            ReplicatesNumber = new ReplicatesNumber(settingValue.ReplicatesNumber, parameterDictionary, userOptionDictionary);
        }

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

        /// <summary>
        /// パラメータファイルに記載する行テキストに変換する。
        /// </summary>
        /// <returns>パラメータ行テキスト</returns>
        public string[] ToParameterFileLines()
        {
            return new[]
            {
                Ploidy.ToParameterFileLine(),
                Parent2PlexNumber.ToParameterFileLine(),
                Bulk1Number.ToParameterFileLine(),
                Bulk2Number.ToParameterFileLine(),
                ReplicatesNumber.ToParameterFileLine()
            };
        }
    }
}
